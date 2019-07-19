class group_rule_input:
    def __init__(self, key_name, output_files,
                 temp_file, group_name):
        '''
        key_name: name to store output under
        output_files: list of lists of files to group together
        temp_files: list of temp files to use to group rules,
                    must have id wildcard
        group_name: name of group to include as
        '''
        self.key_name = key_name
        self.output_files = output_files
        self.temp_file = temp_file
        self.group_name = group_name

    def batch_files(self, batch_size):
        return [self.output_files[i:i+batch_size]
                for i in range(0, len(self.output_files),
                               batch_size)]


class group_rule_output:
    def __init__(self, temp_file, group_name):
        self.temp_file = temp_file
        self.group_name = group_name
        self.files = []

    def get_temp_file(self, i):
        return self.temp_file.format(id=i)

    def get_temp_files(self):
        return [self.get_temp_file(i)
                for i in range(len(self.files))]


def get_batch_files(group_inputs, batch_default, **batches):
    '''
    group_inputs: list of group_rule_input objects
    batch_default: number of files to group together
    batches: keyword arguments of group_name: batch_size to override default
    returns dict of group_rule_outputs keyed on input.key
    '''
    result = {}
    for group_input in group_inputs:
        output = group_rule_output(group_input.temp_file,
                                   group_input.group_name)

        if group_input.key_name in batches:
            output.files = group_input.batch_files(
                batches[group_input.key_name])
        else:
            output.files = group_input.batch_files(batch_default)

        result[group_input.key_name] = output

    return result
