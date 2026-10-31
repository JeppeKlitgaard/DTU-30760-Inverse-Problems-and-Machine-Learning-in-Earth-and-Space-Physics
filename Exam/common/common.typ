#import "exports.typ": *

#let style(doc) = {
  // General setup
  set page(
    columns: 2,
    margin: (
      x: 1.0cm,
      top: 1.5cm,
      bottom: 1.0cm,
    ),
  )

  set text(size: 10pt)

  // Headings
  set heading(numbering: "1.a ")
  show heading.where(level: 1): set text(size: 11pt)
  show heading.where(level: 2): set text(size: 10pt, style: "italic")

  // Prevent line-breaks in the middle of inline math
  show math.equation.where(block: false): it => box(it)

  // Show equation numbers only for referenced equations
  set math.equation(numbering: "(1)")
  show math.equation.where(block: true): it => {
    if it.numbering != none and not it.has("label") {
      math.equation(it.body, numbering: none, block: true)
    } else { it }
  }

  // Show math ref as "(1)"
  show ref: it => {
    let eq = math.equation
    let el = it.element

    if el != none and el.func() == eq {
      link(el.location(), numbering(el.numbering, ..counter(eq).at(el.location())))
    } else {
      it
    }
  }

  // show math.equation.where(block: false): it => {
  //   if it.has("label") and it.at("numbering", default: none) != none {

  //     let cnt = counter(math.equation)
  //     cnt.step()

  //     let supplement = context numbering(it.numbering, ..cnt.get())

  //     math.equation[
  //       #it.body
  //     ]
  //     h(1fr)
  //     supplement

  //   } else { it }
  // }

  // Other formatting
  show: super-T-as-transpose

  doc
}

#let inline-block-math(label, body) = [
  #box[
    #show math.equation: set align(left)
    #show math.equation: set block(spacing: 0pt)
    #math.equation(math.inline(body), block: true)
    #label
  ]
]

#let setup(q-num, q-title, doc) = {
  set document(title: "Exam 30760 — Question " + str(q-num) + ": " + q-title)
  show: style
  doc
}

#let show-page1(doc) = context {
  set page(header: [
    #strong(document.title)
    #h(1fr)
    *Full Notes*
    #line(length: 100%, stroke: 0.7pt)
  ])
  doc
}

#let show-page2(doc) = context {
  set page(header: [
    #strong(document.title)
    #h(1fr)
    *Terse Notes*
    #line(length: 100%)
  ])

  pagebreak(weak: true)
  doc
}

#let unwritten(..args) = [
  #let text = args.pos().at(0)
  #let target = args.pos().at(1, default: none)
  #grid(
    columns: (auto, 1fr, auto),
    align: (left, center, right),
    [
      #set std.text(fill: black.transparentize(25%), style: "italic")
      #text
    ],
    none,
    target
  )
]

#let fill-in-drawing(width, height, caption: none) = [
  #figure(
    box(width: width, height: height, stroke: stroke(paint: black, dash: "dotted")),
    caption: caption,
  )
]