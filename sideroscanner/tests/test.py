from unittest import TestCase
from subprocess import run
from os.path import dirname
import sideroscanner

testfile = f'{dirname(sideroscanner.__file__)}/tests/test.fna'

class TestScan(TestCase):
    def test_is_string(self):
        return run(['sideroscanner', '-i', testfile])