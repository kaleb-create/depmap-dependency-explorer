import importlib
import os
import tempfile
import unittest
from unittest.mock import patch

from stratifier_jobs import StratifierJobs


class StratifierRouteTests(unittest.TestCase):
    def test_saved_build_lifecycle(self):
        with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ, {
            "DB_PATH": os.path.join(directory, "app.db"),
            "DATABASE_URL": "", "DEPMAP_PROVISION_ON_START": "false",
            "AUTH_ENABLED": "false", "OPENAI_API_KEY": "test-only",
        }):
            module = importlib.import_module("app")
            with patch.object(module, "DB_PATH", os.path.join(directory, "app.db")), \
                 patch.object(module, "USING_POSTGRES", False):
                module.init_db()
                jobs = StratifierJobs(module.connect_db, module.db_execute,
                    self.compute, module.save_stratifier, directory)
                with patch.object(module, "stratifier_jobs", jobs):
                    client = module.app.test_client()
                    response = client.post("/stratifiers", json={
                        "prompt": "Test comparison", "request_key": "route-test",
                    })
                    self.assertEqual(response.status_code, 202)
                    job_id = response.json["id"]
                    jobs.thread.join(5)
                    status = client.get(f"/api/stratifier-jobs/{job_id}")
                    self.assertEqual(status.json["status"], "complete")
                    self.assertEqual(client.get("/stratifiers").status_code, 200)
                    self.assertEqual(client.get(status.json["explorer_url"]).status_code, 200)
                    self.assertEqual(client.get("/api/dependency-analysis/custom-1").status_code, 200)
                    analyses = client.get("/api/dependency-summary").json["analyses"]
                    custom = next(item for item in analyses if item["id"] == "custom-1")
                    self.assertNotIn("datasets", custom)
                    self.assertEqual(client.post("/stratifiers", json={"prompt": ""}).status_code, 400)
                    with patch.dict(os.environ, {"OPENAI_API_KEY": ""}):
                        self.assertEqual(client.post("/stratifiers", json={"prompt": "Test"}).status_code, 503)
                    self.assertEqual(client.post("/stratifiers/1/delete").status_code, 302)
                    self.assertEqual(client.get(f"/api/stratifier-jobs/{job_id}").status_code, 404)
                    self.assertEqual(client.get("/api/dependency-analysis/custom-1").status_code, 404)
                    jobs.thread.join(5)

    @staticmethod
    def compute(prompt, progress):
        progress("Testing saved job")
        included = {"positive": ["M0", "M1", "M2"], "negative_n": 3}
        analysis = {
            "label": "Route test", "positive_label": "Positive", "negative_label": "Negative",
            "positive_models": [], "negative_models": [],
            "datasets": {"crispr": [], "rnai": []},
            "included_models": {"crispr": included, "rnai": included},
        }
        source = {
            "dataset": {"name": "Fixture", "provider": "Test", "source_url": ""},
            "retrieval": {"status": "local_metadata"}, "mapping": {},
            "dependency_releases": {"model_metadata": "Test"},
        }
        return analysis, source, {"weaknesses": []}


if __name__ == "__main__":
    unittest.main()
