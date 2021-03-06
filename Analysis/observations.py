import os.path


root = 'photometry/'


def ListClusters():
    """Lists all clusters associated with directories in the project.

    Returns:
        list: A list of strings corresponding to each cluster used in the
              project.

    """
    clusters = []
    if not os.listdir(root):
        print("There are no files in the photometry directory.")
        return clusters

    for cluster in sorted(os.listdir(root)):
        if os.path.isdir(os.path.join(root, cluster)) and \
           cluster[:3] == 'NGC':
            clusters.append(cluster)

    if not clusters:
        print("There are no appropriate cluster directories in the photometry \
               directory.")

    return clusters


def ListDates(cluster):
    """Lists all dates associated with directories in the project for a cluster.

    Args:
        cluster (str): Cluster for which observation dates are found.

    Returns:
        list: A list of strings corresponding to each date used in the
              project for a cluster.

    """
    path = root + cluster + '/'

    dates = []
    if not os.path.isdir(path):
        print("%s is not a valid cluster directory." % path)
        return dates

    if not os.listdir(path):
        print("%s directory does not have any files in it." % cluster)
        return dates

    for date in sorted(os.listdir(path)):
        if os.path.isdir(os.path.join(path, date)) and date[:3] == '201':
            dates.append(date)

    return dates
