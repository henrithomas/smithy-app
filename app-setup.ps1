cd .\smithy
python manage.py makemigrations
python manage.py migrate
python manage.py makemigrations assembly
python manage.py migrate
python manage.py collectstatic
cd ..