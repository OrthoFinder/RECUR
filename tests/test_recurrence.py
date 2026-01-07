def test_recurrence_list(get_recurrence_list):
    expected_recurrence_list = [
        "4 A Q 2 0 A>Q:2 A:5,Q:3",
        "2 K S 2 0 K>S:2 K:6,S:2",
        "6 F Y 2 0 F>Y:2 F:5,Y:3",
        # "10 V I 2 1 V>I:2,I>V:1 V:5,I:2,T:1",
        "10 I V 2 1 I>V:2,V>I:1 V:4,I:3,T:1",
    ]

    assert set(get_recurrence_list) == set(expected_recurrence_list), "RECUR output does not match the expected values"
