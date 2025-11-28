def getBool(dict, key):
    """Return a boolean for 'yes' or 'no' keywords. If the key is not present
    it returns False.

    Args:
        dict (dictionary): The dictionary in which we are searching
        key (key): the keyword we are searching

    Returns:
        Boolean: True if the keywords is set to 'yes'. False otherwise.
    """

    if key in dict:
        boolean = (dict[key] == 'yes')
    else:
        boolean = False
    
    return boolean
