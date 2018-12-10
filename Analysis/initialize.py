def Setup():
    import os.path

    # Set up directories for processing
    setup_directories = [
        'output/',
        'photometry/'
    ]
    for directory in setup_directories:
        if not os.path.exists(directory):
            os.makedirs(directory)
