import pandas as pd
import DataStructures
import os
import sys

csv = ["test.csv",
       "airport.csv",
       "countrycurrency.csv",
       "currencyrates.csv",
       "currencyrates.csv",
       "aircraft.csv"]

missingFiles = []
filesPresent = True

for file in csv:
    if not os.path.isfile(file):
        missingFiles += [file]
        filesPresent = False
        break

if filesPresent:
    test_df = pd.read_csv("test.csv")
    airport_df = pd.read_csv("airport.csv")
    country_currency_df = pd.read_csv("countrycurrency.csv")
    currency_rates_df = pd.read_csv("currencyrates.csv")
    aircraft_df = pd.read_csv("aircraft.csv")

    airport_currencies_df = airport_df.merge(country_currency_df.rename(columns={"name":"Country"}), on="Country")
    merged_df = airport_currencies_df.merge(currency_rates_df.rename(columns={"CurrencyCode":"currency_alphabetic_code"}), on="currency_alphabetic_code")
else:
    print("Error, csv file(s) missing:", missingFiles)
    sys.exit()

g = DataStructures.Graph()

airports = DataStructures.create_airport_list(merged_df)
aircrafts = DataStructures.createAircrafts(aircraft_df)

for i in range(len(test_df)):
    route = DataStructures.getRoute(test_df, i)
    aircraft = test_df["Airplane"][i]
    checkRoute = DataStructures.checkRoute(route, airports)
    if not checkRoute[0]:
        print("Error: Airport,", checkRoute[1], ", doesn't exist in route:", route, "\n")
        continue
    if not DataStructures.checkAircraftExist(aircraft, aircrafts):
        print("Error: Aircraft doesn't exist:", aircraft, "\n")
        continue
    DataStructures.addToAirportAtlas(route, g, airports)
    DataStructures.addAirportDistances(route, g, airports)
    permutations = DataStructures.getPermutations(route)
    print(DataStructures.getCheapestPath(permutations, aircrafts[aircraft].get_maxfuel(), g), "\n")