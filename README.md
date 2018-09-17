# 2048 Neuroevolution

## Introduction
The goal of this project is to create a bot for the game [2048](http://gabrielecirulli.github.io/2048/), using fixed topology neural network evolutionary algorithms, inspired by survival of the fittest. A trial solution, called a *creature* or *individual*, is characterised by its genetic representation, which encodes its behavior. In this project a creature is a neural network, while the genetic representation are the weight matrices of the network. The network takes as input the game's grid and outputs a suitable move. Then a collection of creatures is called a *population*. Neuroevolution can be applied as long as we have a *fitness function* that evaluates each network's performance and *evolution operations*, ways to create creatures based on others. The fitness used is the end game score.

