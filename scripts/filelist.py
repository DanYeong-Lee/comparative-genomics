import os


def mklist(extension, directory='./'):
    file_list = os.listdir(directory)
    target_list = [file.split(extension)[0] for file in file_list if file.endswith(extension)]
    return target_list


def mklist_full(extension, directory='./', absolute=False):
    file_list = os.listdir(directory)
    target_list_full = [file for file in file_list if file.endswith(extension)]
    if absolute is False:
        return target_list_full
    else:
        absolute_target_list = [os.path.abspath(os.path.join(directory, file)) for file in target_list_full]
        return absolute_target_list


