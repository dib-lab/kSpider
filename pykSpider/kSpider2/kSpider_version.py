import os
import requests, json
import sys

# Only update this when releasing stable
MAJOR = 2
MINOR = 1
PATCH = 1

PYPI_PACKAGE = "kSpider"

def is_github_action():
    if "GITHUB_WORKFLOW" in dict(os.environ.items()):
        return True
    else:
        return False    

def get_pypa_dev_latest():
    response = requests.get(f"https://test.pypi.org/pypi/{PYPI_PACKAGE}/json")
    if response.status_code != 404:
        data = json.loads(response.text)
        return data["info"]["version"]
    else:
        # It does not exist yet on test.pypi
        return f"{MAJOR}.{MINOR}.{PATCH}"


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

    # print(f"[DEBUG VERSION] = {version_tag}")
    return version_tag