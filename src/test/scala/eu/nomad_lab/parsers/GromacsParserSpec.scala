package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object GromacsParserTests extends Specification {
  "GromacsParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/aminoacids/md.log", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/aminoacids/md.log", "json") must_== ParseResult.ParseSuccess
    }
  }
}
