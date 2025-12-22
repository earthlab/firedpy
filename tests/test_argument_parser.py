import unittest

from unittest.mock import patch, mock_open
from utilities.argument_parser import FiredpyArgumentParser


class TestFiredpyArgumentParser(unittest.TestCase):

    def setUp(self):
        # Set up a mock params file
        self.mock_params_data = (
            "arg1,Test Prompt,int,5,none\narg2,Do you want to continue?,bool,"
            "y,true|false\narg3,Input a float,float,5.76,none\narg4,Input a "
            "string,str,default,acceptable|string|values\narg5,Input a "
            "password,str,,none"
        )

    def test_load_params(self):
        data = self.mock_params_data
        with patch("builtins.open", mock_open(read_data=data)):
            parser = FiredpyArgumentParser("testfile.txt")
            self.assertIn("arg1", parser.arguments)
            self.assertEqual(parser.arguments["arg1"]["prompt"], "Test Prompt")
            self.assertEqual(parser.arguments["arg1"]["type"], "int")
            self.assertEqual(parser.arguments["arg1"]["last_value"], "5")
            self.assertIsNone(parser.arguments["arg1"]["accepted_values"])

    @patch("builtins.input", return_value="false")
    def test_prompt_for_argument(self, mock_input):
        data = self.mock_params_data
        with patch("builtins.open", mock_open(read_data=data)):
            parser = FiredpyArgumentParser("testfile.txt")
            response = parser.prompt_for_argument("arg2")
            self.assertFalse(response)

    @patch("builtins.input", return_value="")
    def test_prompt_for_argument_default_value(self, mock_input):
        data = self.mock_params_data
        with patch("builtins.open", mock_open(read_data=data)):
            parser = FiredpyArgumentParser("testfile.txt")
            response = parser.prompt_for_argument("arg1")
            self.assertEqual(response, 5)

    @patch("builtins.input", return_value="6.4983")
    def test_prompt_for_float(self, mock_input):
        data = self.mock_params_data
        with patch("builtins.open", mock_open(read_data=data)):
            parser = FiredpyArgumentParser("testfile.txt")
            response = parser.prompt_for_argument("arg3")
            self.assertEqual(response, 6.4983)

    @patch("builtins.input", return_value="acceptable")
    def test_prompt_for_str(self, mock_input):
        data = self.mock_params_data
        with patch("builtins.open", mock_open(read_data=data)):
            parser = FiredpyArgumentParser("testfile.txt")
            response = parser.prompt_for_argument("arg4")
            self.assertEqual(response, "acceptable")

    @patch("builtins.input", side_effect=["not acceptable", "acceptable"])
    def test_prompt_for_str_unacceptable_value(self, mock_input):
        data = self.mock_params_data
        with patch("builtins.open", mock_open(read_data=data)):
            parser = FiredpyArgumentParser("testfile.txt")
            self.assertEqual(parser.arguments["arg4"]["last_value"], "default")
            response = parser.prompt_for_argument("arg4")
            self.assertEqual(2, mock_input.call_count)
            self.assertEqual(response, "acceptable")
            self.assertEqual(
                parser.arguments["arg4"]["last_value"],
                "acceptable"
            )

    @patch("getpass.getpass", return_value="password")
    def test_prompt_for_str_sensitive(self, mock_input):
        data = self.mock_params_data
        with patch("builtins.open", mock_open(read_data=data)):
            parser = FiredpyArgumentParser("testfile.txt")
            self.assertEqual(parser.arguments["arg5"]["last_value"], "")
            response = parser.prompt_for_argument("arg5", sensitive=True)
            self.assertEqual(response, "password")
            self.assertEqual(parser.arguments["arg5"]["last_value"], "")

    def test_prompt_for_non_existent_argument(self):
        data = self.mock_params_data
        with patch("builtins.open", mock_open(read_data=data)):
            parser = FiredpyArgumentParser("testfile.txt")
            with self.assertRaises(ValueError):
                parser.prompt_for_argument("arg_not_exists")


if __name__ == "__main__":
    unittest.main()
