import math
import pandas as pd


# Decorator used to add an equality methods of own objects.
# This equality method checks that all non callable and no dunders attributes are equal
# Moreover 2 attributes that are NaN are also said to be equal


def class_equality_attributes(x):
    def get_non_callable_attributes_and_remove_dunders(obj):
        out = []
        for att in dir(obj):
            if not callable(getattr(obj, att)):
                if len(att) < 4:
                    out.append(att)
                elif (att[:1] != '__') and (att[-2:] != '__'):
                    out.append(att)
        return out

    def iterables_equal(z, y):
        if z == y:
            return True

        if len(z) != len(y):
            return False

        if isinstance(z, list) and isinstance(y, list):
            for i in range(len(z)):
                if z[i] != y[i]:
                    if (isinstance(z[i], list) and isinstance(y[i], list)) or \
                            (isinstance(z[i], dict) and isinstance(y[i], dict)):
                        if not iterables_equal(z[i], y[i]):
                            return False

                    elif isinstance(z[i], float) and isinstance(y[i], float):
                        if not (math.isnan(z[i]) and math.isnan(y[i])):
                            return False
                    else:
                        return False
            return True

        elif isinstance(z, dict) and isinstance(y, dict):
            if z.keys() != y.keys():
                return False
            for i in z.keys():
                if z[i] != y[i]:
                    if (isinstance(z[i], list) and isinstance(y[i], list)) or \
                            (isinstance(z[i], dict) and isinstance(y[i], dict)):
                        if not iterables_equal(z[i], y[i]):
                            return False

                    elif isinstance(z[i], float) and isinstance(y[i], float):
                        if not (math.isnan(z[i]) and math.isnan(y[i])):
                            return False
                    else:
                        return False
            return True

        else:
            print(type(z), type(y))
            return False

    class Wrapper(x):
        def __init__(self, *args, **kwargs):
            x.__init__(self, *args, **kwargs)

        def __eq__(self, other):
            wrap_attributes = get_non_callable_attributes_and_remove_dunders(self)
            other_attributes = get_non_callable_attributes_and_remove_dunders(other)
            # check same number of attributes
            if len(wrap_attributes) != len(other_attributes):
                return False
            # check that they have the same name
            for att in wrap_attributes:
                if att not in other_attributes:
                    return False
            # check that they have the same value
            for att in wrap_attributes:
                if isinstance(getattr(self, att), pd.DataFrame):
                    if not getattr(self, att).equals(getattr(other, att)):
                        return False
                elif getattr(self, att) != getattr(other, att):
                    if isinstance(getattr(self, att), float) and \
                            isinstance(getattr(other, att), float) and \
                            math.isnan(getattr(self, att)) and \
                            math.isnan(getattr(other, att)):
                        continue
                    elif (isinstance(getattr(self, att), list) and isinstance(getattr(other, att), list)) or \
                            (isinstance(getattr(self, att), dict) and isinstance(getattr(other, att), dict)):
                        if iterables_equal(getattr(self, att), getattr(other, att)):
                            continue

                    return False
            return True

    return Wrapper
