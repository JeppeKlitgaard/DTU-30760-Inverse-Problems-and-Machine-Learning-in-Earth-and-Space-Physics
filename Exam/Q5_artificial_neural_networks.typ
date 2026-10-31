#import "common/common.typ": *

#show: setup.with(5, "Artificial Neural Networks")
#show: show-page1

#show heading: set block(above: 0.6em, below: 0.4em)
#show math.equation.where(block: true): set block(above: 0.4em, below: 0.4em)

#let questions = [
= Questions
- What is an *artificial neural network*?
  - Comprised of _neurons_, typically arranged in _layers_, with _weights_ and _biases_ that can be trained to learn a mapping from inputs to outputs.
- How can these be used to *learn a mapping* between example input data and outputs of interest?
- Briefly *give an example application* of their use in Earth and Space Physics
]
#questions

= Setup
- $H$ hidden layers
  - $l ∈ {0, ..., L}, quad L=H+1$
- $M_l$ layer widths
  - $dash(M)_l$ taken to denote inclusion of bias, i.e. index would start at $0$
  - $i ∈ {1, ..., M_l}$ weights
- Input layer: $x_i$, $i∈{1, ..., M_0}$ features
  - $x_0$: bias
- Hidden layers: $z^((j))_k, quad k∈{1, ..., M_j}, j∈{1, ..., H}$
  - $z^((j))_0$: bias for layer $j$
- Output layer: $y_k, quad k∈{1, ..., M_L}$ outputs
- Pre-activations: $a^((j))_k = ∑^(M_(j-1))_(i=1) w^((j))_(k i) z_i + w^((j))_(k 0), j∈{1, .. , L}$
- Activation function: $z_k = h(a_k)$
  - Non-polynomial (thus also non-linear) function (required to satisfy UAT)
  - Examples: Sigmoid, ReLU, Tanh
    // - Sigmoid: $h(a) = 1\/(1 + e^(-a))$
    // - ReLU: $h(a) = max(0, a)$
    // - Tanh: $h(a) = tanh(a)$

= Activation Functions
- Universal Approximation Theorem:
  - Wide: _single hidden layer can approximate any continuous function on compact subset of $ℝ^n$ to arbitrary precision, given enough neurons._
  - Deep: _deep networks can approximate (Lebesque integrable) functions $f : ℝ^n → ℝ$ with $"width" ≤ n+4$._
    - Deep can be exponentially more efficient than shallow

= Example: @fig:ann
- #math.display($y_k (vv(x), vv(w)) = ∑_(j=1)^M w^((2))_(k j) h( ∑_(i=1)^D w^((1))_(j i) x_i + w^((1))_(j 0)) + w^((2))_(k 0)$)
  - Where $vv(w)$ is a _tuple_ of all the parameter matrices $w^((j))_(k i)$
- Can absorb biases by adding $x_0 = z_0 = 1$

= Training
- Training a neural network is just an _inverse problem_, though often formulated as an optimisation problem directly.

== Cost Function
- Define some cost function $E(vv(w))$, expressed as sum of _loss_ functions:
  - $E(vv(w)) = ∑_(n=1)^(N_"data") E_n (vv(w))$
- Often sum of squares for regression:
  - $E(vv(w)) = 1/2 ∑_(n=1)^(N_"data") (y(vv(x)_n, vv(w)) - t_n)^2$
  - For targets $t_n$
- Or cross-entropy for classification (note: not covered):
  - Binary classification (where targets $t_n ∈ {0, 1}$):
    $E(vv(w)) = - ∑_(n=1)^(N_"data") (t_n ln y(vv(x)_n, vv(w)) + (1 - t_n) ln(1 - y(vv(x)_n, vv(w))))$
  - Multiclass classification (for $K$ classes, where $t_(n k)$ is a one-hot encoded target):
    $E(vv(w)) = - ∑_(n=1)^(N_"data") ∑_(k=1)^K t_(n k) ln y_k (vv(x)_n, vv(w))$

== Gradient Descent Methods
- We train by minimising the cost, often using _Stochastic Gradient Descent_ (SGD):
  - $vv(w) ← vv(w) - η ∇ tilde(E)(vv(w))$
  - Where $η$ is the learning rate, $hat(E)$ is loss of a stochastic _mini-batch_.
  - Mini-batching + stochastic advantage:
    - Acts as regularisation, hardens against local minima
    - Pick batch size to fit with parallelisation capabilities

== Backpropagation
- To evaluate gradient efficiently, we use _backpropagation_, which is simply numerical chain rule:
- Change w.r.t weights; $pdv(E_n, w_(i j)) = pdv(E_n, a_j) ⋅ pdv(a_j, w_(i j))$
- Change w.r.t preactivation: $δ_j = pdv(E_n, a_j)$
- $a_j = ∑_i w_(i j) z_i quad ⇒ quad pdv(a_j, w_(i j)) = z_i$
- Chain rule on $δ_j = pdv(E_n, a_j) = ∑_k pdv(E_n, a_k) pdv(a_k, a_j)$
- Plugging in $a_k = sum_j w_(j k) z_j, quad z_j = h(a_j), quad pdv(a_l, a_j) = δ_(l j)$ gives:
$
δ_j
&= sum_k δ_k pdv(, a_j)[∑_l w_(l k) z_l]
&&= sum_k δ_k [∑_l w_(l k) pdv(h(a_l), a_j)]\
&= sum_k δ_k w_(j k) sum_l pdv(a_l, a_j) h'(a_l)
&&= sum_k δ_k w_(j k) h'(a_j)\
&= h'(a_j) sum_k w_(j k) δ_k
$
- Thus we have found _backpropagation formula_: We can find $δ_j$ given by backpropgation of later $δ_k$:
  + Forward propogate to find $a_j, z_j$
  + Evaluate $δ_k$ for outputs using $δ_k = sum_k h'_"out" (a_k) pdv(E_n, y_k)$
  + Backpropagate to find $δ_j$ at hiddne units
  + Evaluate gradient: $pdv(E_n, w_(j i)) = δ_i z_i$

== Regularisation
- Add $L_2$ (ordinary Tikhonov) penality on weights
- Early stopping, similar effect to Tikhonov
- Ensemble methods
- Dropout: Randomly drop units during training. Helps with overfitting.

= Bay. NN: Do Gauss. ≈ on posterior, grad. optim
- Is this actually Laplace Approx?

= CNN
- _Convolve_ input with learnable _filters_
- One _filter_ per _channel_
- Convolution has a _stride_ and input may be _padded_
- _Pooling_ ops downsample (Min, Max, Avg, Norm, ...)
- _Batch norm_ to renormalize batches → stable training
- Often flatten and FCN at end

= Example: Sea Ice Classification
- U-Net architecture, $N=8$:
  - Encoder: N layers, halve spatial dims and double channels
  - Decoder: N layers, double spatial dims and halve channels
  - Skip connections: ≬ enc., dec. layers of same spatial dims
  - ReLU all the way, Softmax at end for classification into 12 classes ($0% — 100%$, and masked)

#let fig-ann-1 = [
  #figure(
    image("ann.png"),
    caption: [
      Single hidden layer ANN.
    ]
  ) <fig:ann>
]

#let fig-ann-2 = [
  #figure(
    image("ann2.png"),
    caption: [
      Two hidden layers.
    ]
  )
]

#v(1fr)
#grid(
  align: bottom,
  columns: 2,
  fig-ann-1,
  fig-ann-2,
)

#show: show-page2
#questions
#LINE
