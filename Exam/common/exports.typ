#import "linalg.typ": mathformatter
#import "@preview/physica:0.9.8": dv, evaluated as eval, pdv, super-T-as-transpose

#let vv = mathformatter(bold: true, upright: false)
#let mm = mathformatter(bold: true, upright: true)
#let wider = h(3em)

#let ddot = math.dot.double
#let dddot = math.dot.triple
#let ul = math.underline
#let Tr = math.op("Tr")
#let Cov = math.op("Cov")

#let LINE = line(length: 100%, stroke: 0.7pt)
#import "@preview/physica:0.9.8": dv, evaluated as eval, pdv, super-T-as-transpose
