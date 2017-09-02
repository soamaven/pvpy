# # content of conftest.py
# from itertools import permutations
#
# def pytest_addoption(parser):
#     parser.addoption("--all", action="store_true",
#         help="run all combinations")
#
# def pytest_generate_tests(metafunc):
#     if 'units' in metafunc.fixturenames:
#         if metafunc.config.getoption('all'):
#             end = 5
#         else:
#             end = 2
#         metafunc.parametrize("x", range(end))
