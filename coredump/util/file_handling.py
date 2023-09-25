import json


def read_json(json_path: str) -> dict:
    """
    Read a JSON file and return its content as a Python dictionary.

    Inputs:
    json_path (str): Path of the JSON file.

    Outputs:
    (dict): Python dictionary containing the content of the JSON file.

    Example:
    Read a JSON file called "config_development.json" located in the
    directory "/home/user/config/" and store the results in a Python
    dictionary called json_data.

    json_data = read_json("home/user/config/config_development.json")
    """

    return json.loads(open(json_path, "r").read())


def read_text_file(text_file_path: str) -> list:
    """
    Read a text file and return its content as a Python list, with each
    element of the list containing one line from the text file.

    Inputs:
    text_file_path (str): Path of the text file.

    Outputs:
    (list): Python list containing the content of the text file, with each
            element of the list containing one line from the text file.

    Example:
    Read a txt file called "requirements.txt" located in the directory
    "/home/user/project_repository/" and store the results in a Python
    list called list_data.

    list_data = read_text_file("home/user/project_repository/requirements.txt")
    """

    return open(text_file_path, "r").read().splitlines()
