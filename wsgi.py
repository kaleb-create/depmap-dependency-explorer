from runtime_data import ensure_depmap_runtime_data


ensure_depmap_runtime_data()

from app import app

if __name__ == "__main__":
    app.run()
