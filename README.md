# tracer

Following along with the University of Utah Ray Tracing course in which, over the semester, the students implement various direct illumination features and, ultimately, implement indirect illumination features.

| Monte Carlo global illumination | Recursive reflection and refraction | Polygonal meshes via Möller-Trumbore |
| --- | --- | --- |
| <img src="https://wilimitis.github.io/assets/img/2018-12-14-project-1.png" width="300"/> | <img src="https://wilimitis.github.io/assets/img/2018-12-10-project-2.png" width="300"/> | <img src="https://wilimitis.github.io/assets/img/2018-12-10-project-1.png" width="300"/> |

Todo
- bdrf experimentation
- acceleration structure experimentation

The above image is a culmination of implementing algorithms to support
- coordinate space transformations: screen, camera, world, object
- recursive scene graph scaling
- naïve lighting: ambient, diffuse, specular
- naïve shadow rays
- recursive reflection and refraction
- polygonal meshes via Möller-Trumbore
- non-adaptive antialiasing
- Monte Carlo global illumination

I immediately realized the difficulties that would result from attempting to follow along _outside_ the class, however didn't give up and worked my way through the first 5 projects before branching off and approaching global illumination. Perhaps the most useful resource I've encounted thus far is Dr. Shirley's [Fundamentals of Computer Graphics](https://www.cs.cornell.edu/~srm/fcg4/).
