import importlib
import pkgutil


def test_import_all_modules():
    # Import all submodules in the resistnet package
    package_name = 'resistnet'
    package = importlib.import_module(package_name)

    for importer, modname, ispkg in pkgutil.walk_packages(
            package.__path__, package_name + '.'):
        importlib.import_module(modname)
