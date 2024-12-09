# smithy
Follow these steps to set up the application on your machine.

1. Install Python 3.11: https://www.python.org/downloads/
2. Install BLAST 2.16: https://blast.ncbi.nlm.nih.gov/Blast.cgi
    - Install and user manual here: https://www.ncbi.nlm.nih.gov/books/NBK569861/
    - Windows install and configuration instructions here: https://www.ncbi.nlm.nih.gov/books/NBK52637/
    - Download the win64 version: ncbi-blast-2.16.0+-win64.exe 
3. Configure BLAST user environment variables to their defaults
    - BLASTDB = "C:\your\db\path\"
    - BLASTDB_LMDB_MAP_SIZE = 1000000
    - Add the BLAST+ bin path to the user's Path environment variable (should be set automatically)
    - e.g. C:\Program Files\NCBI\blast-2.16.0+\bin
4. Copy the provided BLAST (addgene, dnasu, igem) sequence db files into database folder (same as the BLASTDB environment variable)
    - Example db path: "C:\NCBI\db"
5. Add the project path to a PYTHONPATH user environment variable
    - Example: PYTHONPATH = "C:\your\path\smithy-app"
6. Run a pip install of requirements.txt:
    - pip install -r requirements.txt
7. Add the following folders to smithy-app\smithy:
    - media
    - media\images
    - media\csv
8. In command prompt, navigate to smithy-app/smithy
9. Run the app-setup powershell script, or alternatively run steps 10-11.
10. Run the django database migrations:
    - python manage.py makemigrations
    - python manage.py migrate
    - python manage.py makemigrations assembly
    - python manage.py migrate
11. Prepare static files:
    - python manage.py collectstatic
12. Run the app:
    - python manage.py runserver
