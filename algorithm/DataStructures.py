import math
import itertools

# Graph Class
"""
Builds graph representing the costs to travel between airports
"""


class Graph:
    """Methods:
    insertVertex
    getVertex
    insertEdge
    containsVertex
    containsEdge
    getVertices
    getEdge
    getEdgeCount
    getVerticeCount
    lowestCost
    """

    # Graph begins empty with no outs
    # Out will be dictionary of dictionaries so that edges are accessed: out["src"]["dest"]

    def __init__(self):
        self.__out = {}
        self.__edgeCount = 0
        self.__verticeCount = 0

    # Insert a new vertex if it doesn't already exist
    def insertVertex(self, airport=None, conversionRate=None):
        """
        Input: airport and its conversion rate to euro
        Output: New vertex is inserted
        """
        if not self.containsVertex(airport):
            self.__out[airport] = {"conversionRate": conversionRate}
            self.__verticeCount += 1

    # Return a specified vertex
    def getVertex(self, airport):
        """
        Input: An airport name
        Output: An airport vertex
        """
        if self.containsVertex(airport):
            return self.__out[airport]

    # Insert edge if it doesn't already exist
    def insertEdge(self, src, dest, dist=None):
        """
        Input: origin airport, destination airport, distance between
        Output: New edge is inserted between both vertices, accessible from either vertex
        """
        if src != dest and self.containsVertex(src) and self.containsVertex(dest) and not self.containsEdge(src, dest):
            # Insert originally specified edge
            self.__out[src][dest] = dist
            # For going in the opposite direction
            self.__out[dest][src] = self.__out[src][dest]
            self.__edgeCount += 1

    # Check if vertex already exists
    def containsVertex(self, ap):
        """
        Input: An airport
        Output: Return T/F depending on if airport
        """
        return ap in self.getVertices()

    # Check if edge already exists
    def containsEdge(self, src, dest):
        """
        Input: Origin airport, Destination airport
        Ouptut: T/F depending on if edge exists
        """
        return dest in self.__out[src]

    # Return list of vertices
    def getVertices(self):
        """
        Input: None
        Ouput: List of existing vertices
        """
        return list(self.__out.keys())

    # Returns a specific edge between 2 vertices
    def getEdge(self, src, dest):
        """
        Input: Origin airport, Destination airport
        Output: Edge between origin and destination
        """
        if self.containsEdge(src, dest):
            return self.__out[src][dest]

    # Return number of edges
    def getEdgeCount(self):
        """
        Input: None
        Output: Number of edges
        """
        return self.__edgeCount

    # Return number of vertices
    def getVerticeCount(self):
        """
        Input: None
        Output: Number of vertices
        """
        return self.__verticeCount

    def getConversionRate(self, ap):
        """
        Input: Airport
        Output: Conversion rate for airport to euro
        """
        if self.containsVertex(ap):
            return self.__out[ap]["conversionRate"]


# store distance between airport
class hashTable:
    #   create a hash table with empty data and the default size is 0.
    def __init__(self):
        self.__hash_table = [[] for _ in range(5832//4)]
        self.__size = 0
        self.__key = []

    #   inserting key-value pairs into the hash table “O（n）”
    def insert(self, key, value):
        hash_key = hash(key) % len(self.__hash_table)
        key_exists = False
        bucket = self.__hash_table[hash_key]
        #  If the key is already present in the hash table, then update its value with the new one.
        #  Otherwise, inserting new key-value pairs into the hash table, and size plus one
        for i, kv in enumerate(bucket):
            k, v = kv
            if key == k:
                key_exists = True
                break
        if key_exists:
            bucket[i] = ((key, value))
        else:
            bucket.append((key, value))
            self.__key += [key]
            self.__size += 1

    # search if the keys are in the hash table, if it in, return the value   O(n)
    def search(self, key):
        hash_key = hash(key) % len(self.__hash_table)
        bucket = self.__hash_table[hash_key]
        for i, kv in enumerate(bucket):
            k, v = kv
            if key == k:
                return v

    # use to delete keys     O(n)
    def delete(self, key):
        hash_key = hash(key) % len(self.__hash_table)
        key_exists = False
        bucket = self.__hash_table[hash_key]
        # If the key is already in the hash table, then delete that particular key-value pair from the hash table.
        # –Otherwise, no operation is done. and size minus 1
        for i, kv in enumerate(bucket):
            k, v = kv
            if key == k:
                key_exists = True
                break
        if key_exists:
            del bucket[i]
            self.__size -= 1
            self.__key.remove(key)

    # check is the kry and valuw in the hash table  "O(n)"
    def contains(self, key):
        hash_key = hash(key) % len(self.__hash_table)
        bucket = self.__hash_table[hash_key]
        # If the key is already in the hash table, then return true. Otherwise, return false
        for i, kv in enumerate(bucket):
            k, v = kv
            if key == k:
                return True
        return False

    #  O(n^2)
    def keys(self):
        """
        key = []
        for i in self.hash_table:
            if len(i) != 0:
                for j in i:
                    key.append(j[0])
        """
        return self.__key

    # """
    # origianlly, keys use nested for loop to add the keys by optimizing keys, we can now do this in O(1)
    # instesd of O(n^2)
    # """

    # to check id the hash table is empty   O(1)
    def isEmpty(self):
        return self.__size == 0

    # """
    # origianlly, isEmpty use the for loop to check if it's empty,  by optimizing getSize, we can now do this in O(1)
    # instesd of O(n)
    # """

    # O(1)
    def getSize(self):
        return self.__size


# """
# origianlly, getSize use the nested for loop to return the number of keys, we optimized this by store in the
# size directly, increasing it and decrese it  when we need it to , this means the getSize is now O(1) instesd of O(n^2)
# """

class Airport:
    """Methods:
    get_conversion_rate
    get_latitude
    get_longitude
    get_name
    """

    def __init__(self, name, latitude, longitude, conversion_rate=1):
        self.__name = name
        self.__latitude = latitude
        self.__longitude = longitude
        self.__conversion_rate = conversion_rate

    def get_conversion_rate(self):
        return self.__conversion_rate

    def get_latitude(self):
        return self.__latitude

    def get_longitude(self):
        return self.__longitude

    def get_name(self):
        return self.__name

    def get_distance(self, destination):
        # Calculates the distance from a second airport and returns it as a float
        lat1 = self.get_latitude()
        long1 = self.get_longitude()
        lat2 = destination.get_latitude()
        long2 = destination.get_longitude()
        if lat1 == lat2 and long1 == long2:  # the distance between an airport and itself is 0
            return 0
        radius_earth = 6371
        theta1 = long1 * (2 * math.pi) / 360
        theta2 = long2 * (2 * math.pi) / 360
        phi1 = (90 - lat1) * (2 * math.pi) / 360
        phi2 = (90 - lat2) * (2 * math.pi) / 360
        distance = math.acos(
            math.sin(phi1) * math.sin(phi2) * math.cos(abs(theta1 - theta2)) + math.cos(phi1) * math.cos(
                phi2)) * radius_earth
        return math.floor(distance)


# Create class for aircraft

class Aircraft:
    # Methods included in the class:
    """Methods
    __init__:               Object initialisation
    get_maxfuel:            Return the Aircraft's max fuel
    """

    def __init__(self, name='', units='metric', fuelcapacity=0):
        self.__name = name  # Aircraft name
        self.__fuel = 0  # Current fuel in aircraft
        self.__max_Fuel = fuelcapacity  # Max fuel capacity (range is 1km per 1 litre, so aircraft range = max fuel capacity)
        self.__units = units  # units of measurement: metric or imperial

        # If the range is not in metric, convert imperial to metric
        if self.__units != "metric":
            self.__max_Fuel = round((self.__max_Fuel * 1.60934), 2)

    def get_maxfuel(self):
        return self.__max_Fuel

def checkRoute(route, airports):
    fine = True
    notExist = []
    for ap in route:
        if ap not in airports.keys():
            notExist += [ap]
            fine = False
    return fine, notExist

def checkAircraftExist(aircraft, aircrafts):
    return aircraft in aircrafts.keys()

def getRoute(test_df, row):
    """
    Input: dataframe of routes, desired row
    Output: List of airports on the route
    """
    route = list(test_df.drop("Airplane", axis=1).iloc[row])
    return route

def create_airport_list(df):
    airports_list = hashTable() #the data structure to store the airport objects
    rows = df['IATA'].count() #integer for the number of rows in the airport dataframe
    for i in range (0, rows): #creating an airport object from each row
        airports_list.insert(df['IATA'][i], Airport(df['IATA'][i], df['Latitude'][i], df['Longitude'][i], df['toEuro'][i]))
    return airports_list

def createAircrafts(aircraft_df):
    aircrafts = {}
    for i in range(len(aircraft_df)):
        name = aircraft_df["code"][i]
        units = aircraft_df["units"][i]
        fuelcap = aircraft_df["range"][i]
        aircrafts[name] = Aircraft(name, units, fuelcap)
    return aircrafts

def addToAirportAtlas(route, g, airports):
    """
    Input: desired route, Graph g, airports dictionary
    Output: Airports in the route are added to Graph g if not already there
    """
    for ap in route:
        if not g.containsVertex(ap):
            g.insertVertex(ap, airports.search(ap).get_conversion_rate())

def addAirportDistances(route, g, airports):
    """
     Input: desired route, Graph g, airports dictionary
     Output: Edges between airports in the route are added to Graph g if not already there
     """
    for ap in route:
        for ap2 in route:
            if not g.containsEdge(ap, ap2) and not g.containsEdge(ap2, ap):
                g.insertEdge(ap, ap2, airports.search(ap).get_distance(airports.search(ap2)))

# Get permutations where the start airport = the inputted airport, and end airport is = "HOM".
def getPermutations(route):
    desiredPermutations = []
    permutes = list(itertools.permutations(route))
    i = 0
    while permutes[i][0] == route[0]:
        desiredPermutations.append(permutes[i])
        i += 1
    return desiredPermutations


# Calculate the costs of the remaining permutations and return the route with the lowest cost.
def getCheapestPath(permutations, aircraftRange, g):
    # store routes and associated costs in a dictionary
    cost_dict = {}
    start = permutations[0][0]
    for i in permutations:
        i += (start,)
        total = 0
        for j in range(1, len(i)):
            nextJump = g.getEdge(i[j - 1], i[j])
            if nextJump > aircraftRange:
                break

            total += nextJump * g.getConversionRate(i[j-1])
                # Multiply distance by cost

        if i[j] == start:
            cost_dict[total] = i
        else:
            cost_dict[total] = "Impossible"
            
    minimum_cost = min(cost_dict)

    return cost_dict[minimum_cost], round(minimum_cost, 2)
