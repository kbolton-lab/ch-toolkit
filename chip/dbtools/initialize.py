import importlib.resources, sqlite3
from clint.textui import puts, indent

def create_db(db_path, sample_name):
    puts("Constructing database...")
    connection = sqlite3.connect(db_path)
    with connection, indent(4):
        create_tables(connection)
        add_seed_data(connection)
        add_sample_name(connection, sample_name)
    connection.close()

def create_tables(connection):
    puts("Creating tables...")
    source = importlib.resources.files('chip.resources.sql.variantdb').joinpath('schema.sql')
    with importlib.resources.as_file(source) as sql:
        connection.executescript(sql.read_text())

def add_seed_data(connection):
    puts("Populating seed data..")
    resource_dir = importlib.resources.files('chip.resources.sql.variantdb.seed')
    seed_files = [ 'chromosomes.sql',
                   'fp-filter-types.sql',
                   'lofreq-filter-types.sql',
                   'mutect-filter-types.sql',
                   'vardict-filter-types.sql' ]
    with indent(2):
        for f in seed_files:
            puts(f"Seed: {f}")
            source = resource_dir.joinpath(f)
            with importlib.resources.as_file(source) as sql:
                connection.executescript(sql.read_text())

def add_sample_name(connection, sample_name):
    puts(f"Setting sample name to: {sample_name}")
    sql = "INSERT INTO samples (sample_name) VALUES (?)"
    connection.execute(sql, (sample_name,))
