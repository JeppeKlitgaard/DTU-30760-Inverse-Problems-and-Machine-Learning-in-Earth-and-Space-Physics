#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.10": *

#let style(doc) = {
  show: codly-init.with()
  codly(languages: codly-languages)

  show <nonumber>: set heading(numbering: none)
  set math.equation(numbering: "(1)")

  set page(
    numbering: "1/1",
    margin: (
      top: 1.5cm,
      bottom: 1.5cm,
      left: 1.5cm,
      right: 2.5cm,
    )
  )

  doc
}





#set heading(numbering: "1.1: ", supplement: "Part")
#show heading: it => {
  let supp = none
  if it.level == 1 {
    "Part " + counter(heading).display() + it.body
  } else if it.level == 2 {
    "Task " + counter(heading).display() + it.body
  } else if it.level == 3 {
    it.body
  }

}



#set text(
  size: 11pt
)

#let appendix(body) = {
  set heading(numbering: "A.1", supplement: [Appendix])
  counter(heading).update(0)
  body
}
