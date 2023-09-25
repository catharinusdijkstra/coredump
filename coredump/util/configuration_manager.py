from coredump.util.file_handling import (
    read_json,
)


# TODO: Include proper type hinting and documentation here.
class ConfigurationManager:
    def __init__(self, configuration_file: str) -> None:
        self.configuration_file = configuration_file

    def get_configuration(self):
        self.configuration = read_json(self.configuration_file)
        return self.configuration
