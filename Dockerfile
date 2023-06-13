FROM harm:latest

COPY ./launch_harm.sh /code/launch_harm.sh
RUN chmod 777 /code/launch_harm.sh
WORKDIR /code

# This is the command that will run your model
ENTRYPOINT ["/code/launch_harm.sh"]


