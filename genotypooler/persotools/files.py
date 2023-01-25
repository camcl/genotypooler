import os
import typing
import re
import pathlib


'''
Personal tools for files management
'''


class FilePathPattern(str):
    """
    Duck class for experimentation.
    """
    def __new__(cls, fpath):
        pattern = re.compile(r'^.*/(?!.*/)(.*)$')
        if pattern.match(fpath) is None:
            fpath = None
        return super().__new__(cls, fpath)


FilePath = typing.NewType('FilePath', typing.Union[pathlib.Path, os.PathLike])


def delete_file(file_path: FilePath) -> None:
    """
    Deletes an existing local file
    :param file_path: string
    :return: -
    """
    if os.path.exists(file_path):
        os.remove(file_path)
    else:
        print("{}: the file does not exists".format(file_path))


def file_size(fname: str) -> float:
    """

    :param fname:
    :return:
    """
    statinfo = os.stat(fname)
    return statinfo.st_size


def dir_size(dirpath: typing.Union[str, pathlib.Path]) -> float:
    node = os.path.split(dirpath)[0]
    os.chdir(node)
    parent = os.path.split(dirpath)[1]
    sz = 0
    for walker in os.walk(parent):
        for file in walker[2]:
            try:
                sz =+ os.stat(os.path.join(node, walker[0], file)).st_size
            except:
                pass
    return sz


def mkdir(dir_path: FilePath) -> None:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def check_file_creation(wd: FilePath, fname: FilePath) -> bool:
    """

    :param wd:
    :param fname:
    :return:
    """
    file_path = os.path.join(wd, fname)
    return os.path.exists(file_path)


def locate(func):
    def function_wrapper(*args, **kwargs):
        print('Current working directory: ', os.getcwd())
        func(*args, **kwargs)
    return function_wrapper
