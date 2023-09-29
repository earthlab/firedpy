from typing import List, Any
import getpass


class FiredpyArgumentParser:

    def __init__(self, params_file):
        self.params_file = params_file
        self.arguments = self._load_params()

    def _load_params(self):
        args = {}
        try:
            with open(self.params_file, 'r') as file:
                for line in file:
                    name, prompt, arg_type, last_value, accepted_values = line.strip().split(',')
                    args[name] = {
                        'prompt': prompt.replace('\\n', '\n'),  # Convert any literals to newlines
                        'type': arg_type,
                        'last_value': last_value if last_value != 'none' else None,
                        'accepted_values': None if accepted_values == 'none' else accepted_values.split('|')
                    }
        except FileNotFoundError:
            pass

        return args

    def _save_params(self):
        # TODO: Think about making this safer so any errors don't overwrite the only existing file
        with open(self.params_file, 'w') as file:
            for name, data in self.arguments.items():
                # Convert newlines back to literals before writing
                prompt_literal = data['prompt'].replace('\n', '\\n')
                accepted_values = 'none' if data['accepted_values'] is None else data['accepted_values']
                file.write(f"{name},{prompt_literal},{data['type']},{data['last_value']},{accepted_values}\n")
        self._load_params()

    def prompt_for_argument(self, arg_name, prompt_override: str = None, accepted_value_override: List[Any] = None,
                            sensitive: bool = False):
        if arg_name not in self.arguments:
            raise ValueError(f"Argument {arg_name} not found in params file.")

        arg_data = self.arguments[arg_name]
        last_value = arg_data['last_value']
        prompt = arg_data['prompt'] if prompt_override is None else prompt_override
        accepted_values = arg_data['accepted_values'] if accepted_value_override is None else accepted_value_override

        if sensitive:
            user_input = getpass.getpass(prompt)
        else:
            user_input = input(f"{prompt} [Default: {last_value}]: ")

        # If no input is given, use the last value.
        if not user_input:
            user_input = last_value

        if accepted_values is not None and user_input not in accepted_values:
            print(f'{user_input} is not an acceptable value. Acceptable values are:\n {accepted_values}')
            return self.prompt_for_argument(arg_name, prompt_override)

        # Convert the input to the appropriate type.
        if arg_data['type'] == 'int':
            user_input = int(user_input)
        elif arg_data['type'] == 'float':
            user_input = float(user_input)
        elif arg_data['type'] == 'bool':
            user_input = user_input.lower() in ['y', 'true', 'yes']

        # Update the last used value.
        if not sensitive:
            self.arguments[arg_name]['last_value'] = user_input

        # Save updated parameters to file.
        self._save_params()

        return user_input

