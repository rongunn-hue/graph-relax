import tempfile,unittest
from pathlib import Path
import render_iamp_common_metric_direct_clean as r
class Tests(unittest.TestCase):
 def test_render_only_and_annotation_count(self):
  source=Path("output/graph_first");files,maps=r.run(source);self.assertEqual(len(files),5);self.assertEqual(len(maps["start"]),34);self.assertEqual(len(maps["final"]),34)
  for p in files[:4]:
   t=p.read_text();self.assertNotIn("d/D=",t);self.assertNotIn("EVENT TABLE",t)
  t=files[-1].read_text();self.assertIn("EVENT TABLE",t);self.assertIn("TPO: own TPO–EOUT excluded",t);self.assertIn("TPB: NODE_B–TPB",t)
if __name__=="__main__":unittest.main()
