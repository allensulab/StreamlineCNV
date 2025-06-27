# Building the image yourself

To create your own Docker image using these files, [follow the instructions to install Docker](https://docs.docker.com/desktop/).
After installing Docker, navigate to this folder and run the following:

```
docker build -t example/tag -f Dockerfile  
```

This will create an image with the tag you provided and store it locally on your machine.
