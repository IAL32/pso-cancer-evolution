# Particle Swarm Optimization for Cancer Evolution

## Pseuco-Code

```
for each particle
    initialize position particle
    initialize velocity particle
end for

// iteration
k = 1
do
    for each particle
        calculate fitness value
        if < the fitness value is better than the best fitness value (pbest) in history >
            set current value as the new pbest
        end if
    end for

    choose the particle with the best fitness value of all the particles as the gbest

    for each particle
        calculate particle velocity according to:
            v_i(k + 1) = w * v_i(k) + c1 * rand(p_i - x_i) + c2 * rand(p_g - x_i)
        calculate particle position according to:
            x_i(k + 1) = x_i(k) + v_i(k + 1)
    end for
    k = k + 1
while maximum interations or minimum error criteria is not attained
```

## The General Idea

The aim is to find an optimal tree given the input matrix.

## Setup

The setup is very easy, just install every package included in the project and run with the required arguments!

I recommend creating a `virtualenv` with `python3`, because that's what I used:
```shell
$ virtualenv env -p `which python 3`
$ source ./env/bin/activate
```

Now you should be in a situation like this:

```shell
(env) $ _ 
```

And then install every package necessary:

```shell
(env) $ pip install graphviz ete3
```

Done!

## Examples

Here's a list of examples that you can try:

```shell
(env) $ python main.py --infile "data/easy/3.txt" --mutations 8 --iterations 10 --mutfile "data/easy/3_mut.txt"
```

```shell
(env) $ python main.py --infile "data/scg_gawad/pat1.txt" --mutations 20 --mutfile "data/scg_gawad/pat1_mut.txt" --particles 10 --iterations 500
```

# Diestel - Graph Theory
# Korte Vygen - Combinatorial Optimization