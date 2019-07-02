import unittest
from get_triplet_repeats_3_alleles import *


class test_triplet_repeats(unittest.TestCase):


    def test_get_triplets_table (self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        self.assertEqual(len(triplets_table),7)
        self.assertEqual(len(triplets),12)
        self.assertEqual(triplets_table.iloc[0,3],406.0)        
        self.assertEqual(triplets_table.iloc[4,4],258.9)


    #def test_match_control_samples_with_references(self):
        triplets,triplets_table=get_triplets_table("gene1","tester1")
        controls,continue_program=match_control_samples_with_references(triplets, "gene1")
        self.assertEqual(continue_program, "yes")
        self.assertEqual(len(controls),4)
        self.assertEqual(controls.iloc[0,4],'51')
        self.assertEqual(controls.iloc[3,4],'29')

    #def test_find_closest_control_peak_to_sample_peaks(self):



    #def test_get_number_of_triplet_repeats(self):



    #def test_format_columns(self):







if __name__=='__main__':
    unittest.main()
