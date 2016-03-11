import pickle

tuple_dict = pickle.load( open( "dicts_baby.p", "rb" ))

transdom_intra_distances = tuple_dict[0]
transdom_inter_distances = tuple_dict[1]
transdom_aa_freq = tuple_dict[2]


for key,value in sorted(transdom_inter_distances.items()):
    print(key)
    for key, value in sorted(value.items()):
        string = key ,value
        print(string, end = "\t")
