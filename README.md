# dystopic-277

## Julia
[I suggest taking a look at the style guide.](https://docs.julialang.org/en/v1/manual/style-guide/index.html) Style isn't just about making things pretty (although that should be reason enough!) it's also about making things clear and making communication easy.

In addition, I'd emphasize the following:
 - Write functions. Don't write scripts. Scipts are bad. 
   - Scripts are not reusable, meaning more work later on when you make changes
   - Scripts are not clear, a script is long an unruly, functions are short and sweet
   - Global behavior; never a good thing in Julia. Julia gets its performance from the concept of local scope, which functions provide. You can read about type stability if you're interested in more about performance, otherwise I'll just take care of it.
   - functions have names that say what they're doing
   - functions can be documented easily.
   - functions place scope boundaries between different conceptual tasks
   - functions improve debugability.
   
So anyway you should write functions to do small tasks and then put them together to do bigger things (true in every language!!). It will make our code modular and assemblable, which is something we'll need in order to do a variety of different tests. Take a look at what I did with the code in graph_environment.jl for ideas on how. I'm sure it could be even shorter and cleaner but whatever. 


## Main Aspects (the two are mostly independent)

 - __Patrol__: scan a region as best as possible
   - Define successful detection
 - __Pursuit__: kill kill kill. No survivors(leave one survivor to tell their friends of the horrors)
   - Define successful capture
   - Define failure condition
 
 ## Necessary tools:
 
  - Environment representation (options) 
      1) A graph with underlying geometry (Ethan uploaded file for this - working on tool to plot geometry from it tomorrow)
      2) Geometry with an underlying graph
      3) __A graph with no associated geometry and some extra plotting machinery__ -- looks like we're going with this option -- I would say it's more like #1 --jv
    
 - Graph traversal
   - Dynamics

 - Agent Controls
   - Stochasticity
   - MDPs? POMDPs?
   - Communication/distribution
 
 - Intruder Controls
   - Omnicient? 
   - Relative properties?
   - Multiple intruders?
   - Intrusion probability 
   - Goal
   
