# tracer

A multi-threaded, batch renderer built on top of the University of Utah Ray Tracing course sample code.
Primary reference: Dr. Shirley's [Fundamentals of Computer Graphics](https://www.cs.cornell.edu/~srm/fcg4/).

| Monte Carlo global illumination | Recursive reflection and refraction | Polygonal meshes via Möller-Trumbore |
| --- | --- | --- |
| <img src="https://wilimitis.github.io/assets/img/2018-12-14-project-1.png" width="300"/> | <img src="https://wilimitis.github.io/assets/img/2018-12-10-project-2.png" width="300"/> | <img src="https://wilimitis.github.io/assets/img/2018-12-10-project-1.png" width="300"/> |

Todo
- bdrf experimentation
- acceleration structure experimentation

The above image is a culmination of implementing algorithms to support
- coordinate space transformations: screen, camera, world, object
- recursive scene graph scaling
- naïve multi-threading
- naïve lighting: ambient, diffuse, specular
- naïve shadow rays
- non-adaptive antialiasing
- polygonal meshes via Möller-Trumbore
- recursive reflection and refraction
- Monte Carlo global illumination
