
class Msg:

    @staticmethod
    def msg_data_gen(func: str, data_missing: str, func_call: list) -> str:
        """
        Prints the error string for missing data (i.e. not yet generated and the list of functions to generate this data)
        Parameters
        ----------
        data_missing
        func_call

        Returns
        -------

        """
        func_str = '\n'.join(func_call)
        return f'ERR in function: {func}\n' \
               f'You have not yet generated: {data_missing}. Please look at the following functions to generate that ' \
               f'data: \n{func_str}'

    @staticmethod
    def msg_data_type(func: str, data: object, data_type_req: str) -> str:
        return f'ERR in function: {func}\n' \
               f'Your data was of type: {type(data)}, variable passed: {data}.\n' \
               f'We require data to be of type: {data_type_req}'

    @staticmethod
    def msg_arg_err(func: str, arg_name: str, value: object, req: list) -> str:
        req_str = '\n'.join(req)
        return f'ERR in function: {func}\n' \
               f'For argument ({arg_name}) you passed: {value}, while {arg_name} must be one of: \n{req_str}.'

    @staticmethod
    def msg_file_not_found(func: str, filename: str) -> str:
        return f'ERR in function: {func}\n' \
               f'File: {filename} does not exist.'