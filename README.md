# smithy
Follow these steps to set up the application on your machine.

1. Install Python 3.8: https://www.python.org/downloads/
2. Install BLAST 2.11 (newer versions may work): https://blast.ncbi.nlm.nih.gov/Blast.cgi
    - Install and user manual here: https://www.ncbi.nlm.nih.gov/books/NBK569861/
    - Windows install and configuration instructions here: https://www.ncbi.nlm.nih.gov/books/NBK52637/
3. Configure BLAST user environment variables to their defaults
    - BLASTDB = C:\your\db\path\
    - BLASTDB_LMDB_MAP_SIZE = 1000000
    - Add the BLAST+ bin path to the user's Path environment variable
4. Copy the provided BLAST (addgene, dnasu, igem) sequence db files into database folder (same as the BLASTDB environment variable)
    - Example db path: C:\NCBI\db
5. Add the project path to a PYTHONPATH user environment variable
    - PYTHONPATH = C:\...\smithy-app
6. Run a pip install of requirements.txt:
    - pip install -r requirements.txt
7. In command prompt, navigate to smithy-app/smithy
8. Run the django database migrations:
    - python manage.py makemigrations
    - python manage.py migrate
    - python manage.py makemigrations assembly
    - python manage.py migrate
9. Run the app:
    - python manage.py runserver
