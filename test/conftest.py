"""
A home for common test fixtures
"""

import os
import pytest


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')

# file paths for additional findings
ACMG_1 = os.path.join(INPUT, 'acmg_test_1.csv')
ACMG_1_EXP = os.path.join(INPUT, 'acmg_1_expected.json')
ACMG_2 = os.path.join(INPUT, 'acmg_test_2.csv')
ACMG_2_EXP = os.path.join(INPUT, 'acmg_2_expected.json')
ACMG_3 = os.path.join(INPUT, 'acmg_test_3.csv')
ACMG_3_EXP = os.path.join(INPUT, 'acmg_3_expected.json')
ACMG_4 = os.path.join(INPUT, 'acmg_test_4.csv')
ACMG_4_EXP = os.path.join(INPUT, 'acmg_4_expected.json')


@pytest.fixture(name='acmg_1', scope='session')
def fixture_acmg_1() -> str:
    """
    :return: STR path to the ACMG file
    """
    return ACMG_1


@pytest.fixture(name='acmg_1_exp', scope='session')
def fixture_acmg_1_expected() -> str:
    """
    :return: STR path to the ACMG file
    """
    return ACMG_1_EXP


@pytest.fixture(name='acmg_2', scope='session')
def fixture_acmg_2() -> str:
    """
    :return: STR path to the ACMG file
    """
    return ACMG_2


@pytest.fixture(name='acmg_2_exp', scope='session')
def fixture_acmg_2_expected() -> str:
    """
    :return: STR path to the ACMG file
    """
    return ACMG_2_EXP


@pytest.fixture(name='acmg_3', scope='session')
def fixture_acmg_3() -> str:
    """
    :return: STR path to the ACMG file
    """
    return ACMG_3


@pytest.fixture(name='acmg_3_exp', scope='session')
def fixture_acmg_3_expected() -> str:
    """
    :return: STR path to the ACMG file
    """
    return ACMG_3_EXP


@pytest.fixture(name='acmg_4', scope='session')
def fixture_acmg_4() -> str:
    """
    :return: STR path to the ACMG file
    """
    return ACMG_4


@pytest.fixture(name='acmg_4_exp', scope='session')
def fixture_acmg_4_expected() -> str:
    """
    :return: STR path to the ACMG file
    """
    return ACMG_4_EXP
