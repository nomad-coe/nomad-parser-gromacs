package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object GromacsParserTests extends Specification {
  "GromacsParserTest" >> {
//    "[aminoacids] test with json-events" >> {
//      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/aminoacids/md.log", "json-events") must_== ParseResult.ParseSuccess
//    }
//    "[aminoacids] test with json" >> {
//      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/aminoacids/md.log", "json") must_== ParseResult.ParseSuccess
//    }
//    "[argon] test with json-events" >> {
//      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/argon/md.log", "json-events") must_== ParseResult.ParseSuccess
//    }
//    "[argon] test with json" >> {
//      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/argon/md.log", "json") must_== ParseResult.ParseSuccess
//    }
//    "[water] test with json-events" >> {
//      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/water/md.log", "json-events") must_== ParseResult.ParseSuccess
//    }
//    "[water] test with json" >> {
//      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/water/md.log", "json") must_== ParseResult.ParseSuccess
//    }
//    "[tip4p] test with json-events" >> {
//      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/tip4p/md.log", "json-events") must_== ParseResult.ParseSuccess
//    }
//    "[tip4p] test with json" >> {
//      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/tip4p/md.log", "json") must_== ParseResult.ParseSuccess
//    }
    "[Fe] test with json-events" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/fe_test/md.log", "json-events") must_== ParseResult.ParseSuccess
    }
    "[Fe] test with json" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/fe_test/md.log", "json") must_== ParseResult.ParseSuccess
    }
  }
}
