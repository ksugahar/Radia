import unittest
from phase_memory import instrument


class PhaseMemoryTest(unittest.TestCase):
    def test_success_preserves_result_and_arguments(self):
        events=[]
        result=object()
        def solve(x, *, y):
            self.assertEqual((x,y),(3,4))
            return result
        fn=instrument(solve,lambda *a,**k:events.append(k),lambda:123,interval=1)
        self.assertIs(fn(3,y=4),result)
        self.assertTrue(events[0]['completed'])
        self.assertEqual(events[0]['sampled_max_process_rss_bytes'],123)
        self.assertEqual(events[0]['samples'],2)

    def test_failure_is_not_pass(self):
        events=[]
        def solve(): raise ValueError('solver failed')
        with self.assertRaisesRegex(ValueError,'solver failed'):
            instrument(solve,lambda *a,**k:events.append(k),lambda:123)()
        self.assertFalse(events[0]['completed'])

    def test_measurement_failure_is_explicit(self):
        events=[]
        def rss(): raise RuntimeError('unavailable')
        self.assertEqual(instrument(lambda:9,lambda *a,**k:events.append(k),rss)(),9)
        self.assertIsNone(events[0]['sampled_max_process_rss_bytes'])
        self.assertTrue(events[0]['errors'])


if __name__=='__main__': unittest.main()
