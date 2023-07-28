# Inside the ch-toolkit repository

        # do some work...
        git commit             # commit things into the repository
        git tag -a v2.2.1      # make a tag on that commit v<version-number>
        git push origin        # push all the commits to github
        git push origin --tags # push all the tags to github
Now if you go to https://github.com/kbolton-lab/ch-toolkit/tags you should see the tag just pushed

## [Semantic Versioning 2.0.0](https://semver.org/)

Given a version number MAJOR.MINOR.PATCH, increment the:
1. MAJOR version when you make incompatible API changes
2. MINOR version when you add functionality in a backward compatible manner
3. PATCH version when you make backward compatible bug fixes

# Inside the docker-toolkit repository
    
        cd docker-toolkit/           # After cloning repo: https://github.com/kbolton-lab/docker-toolkit
        cd ch-toolkit/
        vim -p Dockerfile README.md  # change things to the latest "tag" version
        docker images
        docker build -t ch-toolkit:v2.2.1 .
        docker run -i -t --rm ch-toolkit:v2.2.1 /bin/bash
        docker images
        docker tag d53f12a9a39b indraniel/ch-toolkit:v2.2.1
        docker images
        docker push indraniel/ch-toolkit:v2.2.1