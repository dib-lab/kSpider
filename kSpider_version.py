import os
import json
import urllib.request
import sys

#TODO there are two files with kSpider_version, remove one.

# Only update this when releasing stable
MAJOR = 2
MINOR = 3
PATCH = 2

PYPI_PACKAGE = "kSpider"

def is_github_action():
    if "GITHUB_WORKFLOW" in dict(os.environ.items()):
        return True
    else:
        return False    

def get_pypa_dev_latest():
    # Should handle if the package does not exist.
    with urllib.request.urlopen("https://test.pypi.org/pypi/kSpider/json") as url:
        data = json.loads(url.read().decode())
        return data["info"]["version"]

def increment_patch_version(patch_version):
    patch_version = patch_version.split('.')
    modified_version = f"{patch_version[0]}.{patch_version[1]}.{int(patch_version[2]) + 1}"
    if len(patch_version) == 4:
        modified_version += f".{patch_version[3]}"
    else:
        modified_version += f".dev0"
    return modified_version


def get_version():
    
    dev_version = f"{MAJOR}.{MINOR}.{PATCH}.dev0"
    release_version = f"{MAJOR}.{MINOR}.{PATCH}"
    
    version_tag = str()

    if "BRANCH_NAME" in os.environ:
        print(f"BRANCH_NAME: {os.environ['BRANCH_NAME']}")

        if "master" in os.environ["BRANCH_NAME"]:
            version_tag = release_version
        else:
            test_pypa_latest_version = get_pypa_dev_latest()
            version_tag = increment_patch_version(test_pypa_latest_version)

        # Running on local machine
    else:
        version_tag = dev_version

    print(f"[DEBUG VERSION] = {version_tag}")
    return version_tag
