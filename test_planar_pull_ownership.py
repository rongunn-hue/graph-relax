import unittest
import planar_pull_ownership as ownership

class PullOwnershipTests(unittest.TestCase):
    def setUp(self):
        # The same degree-2 vertex is boundary-incident and FREE-FOR-FACE on
        # both sides of a two-edge cycle boundary representation.
        self.walks={"face:010":["a","v","b"],"face:020":["b","v","a"]}
        self.edges=[("a","v"),("v","b"),("a","b")]

    def test_shared_incidence_exclusive_lowest_layer_ownership(self):
        boundary,free,_=ownership.derive_face_incidence_and_free_for_face(self.walks,self.edges)
        self.assertEqual(["face:010","face:020"],boundary["v"]);self.assertEqual(["face:010","face:020"],free["v"])
        owners,cells,_=ownership.assign_pull_ownership(boundary,free,{"face:010":2,"face:020":1})
        self.assertEqual("face:020",owners["v"]);self.assertEqual(1,sum("v" in members for members in cells.values()))

    def test_same_layer_lexical_tie_and_single_destination(self):
        boundary,free,_=ownership.derive_face_incidence_and_free_for_face(self.walks,self.edges)
        owners,cells,_=ownership.assign_pull_ownership(boundary,free,{"face:010":1,"face:020":1})
        self.assertEqual("face:010",owners["v"])
        coords={"a":(0.,0.),"v":(1.,0.),"b":(2.,0.)};dest,factors=ownership.calculate_owned_destinations(coords,(0.,1.),owners,{"face:010":.9,"face:020":.8})
        self.assertEqual(set(owners),set(dest));self.assertEqual(.9,factors["v"])

    def test_plant_shared_vertex_once_and_reconnect_all_edges(self):
        coords={"a":(0.,0.),"v":(1.,0.),"b":(2.,0.)};new=ownership.plant_destinations(coords,{"v":(1.,1.)})
        segments=ownership.rebuild_edge_segments(new,self.edges);ownership.assert_reconnected(segments,new,self.edges)
        self.assertEqual([1.,1.],segments[0]["end"]);self.assertEqual([1.,1.],segments[1]["start"]);self.assertEqual([0.,0.],segments[2]["start"])

if __name__=='__main__':unittest.main()
