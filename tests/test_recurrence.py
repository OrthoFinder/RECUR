import recur

def test_recurrence_list(get_recurrence_list):
    expected_recurrence_list = ["2 K S 2 0 K>S:2 K:6,S:2",
                                "5 R Q 2 0 R>Q:2 Q:5,R:3",
                                "6 F Y 2 0 F>Y:2 F:5,Y:3"]
    assert all(rec_list in expected_recurrence_list for rec_list in get_recurrence_list), "RecFinder output does not match the expected values"
