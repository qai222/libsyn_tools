import os

# paths
UTILS_DIR = os.path.dirname(os.path.abspath(__file__))

# askcos url
ASKCOS_URL_TXT = os.path.join(UTILS_DIR, "askcos_url.txt")
with open(ASKCOS_URL_TXT, "r") as f:
    ASKCOS_URL = f.read().strip()
