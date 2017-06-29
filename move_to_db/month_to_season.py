def get_season_by_month(month):

    month_to_season = {"winter" : [1, 2, 11, 12],
                       "summer" : [5, 6, 7, 8],
                       "equinox" : [3, 4, 9, 10]}

    for ky in month_to_season.keys():
        if month in month_to_season[ky]:
            return ky
