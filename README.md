# smithy
Follow these steps to set up the application on your machine.

1. Install Python 3.8: https://www.python.org/downloads/
2. Install BLAST 2.11 (newer versions may work): https://blast.ncbi.nlm.nih.gov/Blast.cgi
3. Configure BLAST environment variables to their defaults
4. Copy BLAST sequence db files into database folder (same as the BLASTDB environment variable)
6. Run a pip install of requirements.txt
7. In command prompt, navigate to smithy-app/smithy 
8. Run the django database migrations:
    - python manage.py makemigrations
    - python manage.py migrate
    - python manage.py makemigrations assembly
    - python manage.py migrate
9. Run the app:
    - python manage.py runserver
