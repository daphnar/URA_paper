def m2inch(value,unit='mm'):
    if unit=='cm':
        return value/2.54
    if unit=='mm':
        return value/10/2.54
def one_column():
    return m2inch(89)

def two_columns():
    return m2inch(183)

def full_page():
    return m2inch(247)