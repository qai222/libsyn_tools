from __future__ import annotations

import types
import typing
from typing import Type, Optional, get_origin

from loguru import logger
from pydantic import create_model, BaseModel
from pydantic.fields import FieldInfo
from pydantic_core import PydanticUndefined
from twa.data_model.base_ontology import BaseClass
from twa.data_model.base_ontology import DatatypeProperty
from twa.data_model.utils import construct_rdf_type


def pydantic_to_ontology_class(pydantic_class: Type[BaseModel], base_ontology_class: Type[BaseClass]):
    """
    create an ontology class from a pydantic class
    # TODO right now this only considers data properties
    # TODO would it be good to have data property names as "has_<field name>"?

    originally by Jiaru

    caveats:
    - it assumes the data property is defined under the same namespace of the concept, which might not be true if someone defined properties with the same name in different ontologies, i.e.
    ```
    class MyClass1(BaseClass):
        rdfs_isDefinedBy = MyOntology1
        myData: MyData[str]

    class MyClass2(BaseClass):
        rdfs_isDefinedBy = MyOntology1
        myData: myOntology2.MyData[int] # if MyData here is intended to be defined in another ontology other than MyOntology1
    ```
    - probably doesn't work for fields defined with more than two possible type hints, i.e.
    ```
        myField: str | float | None
    ```
    - the code assumes the fields in the pydantic class (data properties in the output ontology class) are defined with camelCase or snake_case (the first letter is always lower letter)

    :param pydantic_class:
    :param base_ontology_class:
    :return:
    """

    field_definitions = dict()
    for k, v in pydantic_class.model_fields.items():
        logger.info(f"working on pydantic field {k}")
        if k == "identifier":
            continue
        if k in ("rdfs_isDefinedBy", "instance_iri"):
            raise RuntimeError(f"a key of the pydantic class is {k}, this would likely cause errors for ontology class")
        # check if data property already exist
        data_property_predicate_iri = construct_rdf_type(
            base_ontology_class.rdfs_isDefinedBy.namespace_iri,
            k[:1].lower() + k[1:]
        )

        # this is ugly but somehow `data_property_lookup` could be None instead of `{}`
        if base_ontology_class.rdfs_isDefinedBy.data_property_lookup is None:
            data_property_exists = False
        elif data_property_predicate_iri in base_ontology_class.rdfs_isDefinedBy.data_property_lookup:
            data_property_exists = True
        else:
            data_property_exists = False

        if data_property_exists:
            data_property_clz = base_ontology_class.rdfs_isDefinedBy.data_property_lookup[data_property_predicate_iri]
        else:
            data_property_clz = DatatypeProperty.create_from_base(k[:1].upper() + k[1:],
                                                                  base_ontology_class.rdfs_isDefinedBy)

        # parse annotation of the field v
        v: FieldInfo
        if isinstance(v.annotation, (typing._UnionGenericAlias, types.UnionType)):
            type_args = set(v.annotation.__args__)
            if types.NoneType in type_args:
                type_args.remove(types.NoneType)
                range_type_lst = list(type_args)
                if len(range_type_lst) > 1:
                    raise NotImplementedError(
                        'More than one possible range for datatype property is not yet supported.')
                data_property_type = Optional[data_property_clz[range_type_lst[0]]]
            else:
                raise NotImplementedError('More than one possible range for datatype property is not yet supported.')
        else:
            data_property_type = data_property_clz[v.annotation]

        # TODO when `data_property_type` is a derivative of list, `create_model` errors out
        #  twa\data_model\base_ontology.py", line 644, in __get_pydantic_core_schema__
        #      if issubclass(tp, BaseClass):
        #         ^^^^^^^^^^^^^^^^^^^^^^^^^
        #    File "<frozen abc>", line 123, in __subclasscheck__
        #  TypeError: issubclass() arg 1 must be a class
        #  the following escape such fields
        if get_origin(v.annotation) in (list, tuple):
            logger.critical(
                f"in creating ontology class for {pydantic_class.__name__} "
                f"we are skipping field: {k} with range annotation: {v.annotation}"
            )
            continue
        logger.info(f'Complete field annotation will be: {data_property_type}\nDefault is {v.default}')
        field_definitions[k] = (data_property_type, v.default if v.default is not PydanticUndefined else ...)

    ontology_class = create_model(
        pydantic_class.__name__,
        __base__=base_ontology_class,
        **field_definitions
    )

    return ontology_class
