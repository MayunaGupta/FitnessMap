import sqlite3
import pandas as pd
conn = sqlite3.connect("feba.db")
df = pd.read_sql("SELECT * FROM fitness LIMIT 5", conn)
