import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ngl_resum",
    version="0.0.1",
    author="Marcel Balsiger",
    author_email="marcel.balsiger@hotmail.com",
    description="Resummation of non-global logarithms at leading logarithmic accuracy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MarcelBalsiger/ngl_resum",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
