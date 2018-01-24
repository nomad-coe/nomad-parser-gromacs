package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object GromacsParserTests extends Specification {
  "GromacsParserTest_1" >> {
    "[aminoacids] test with json-events" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/aminoacids/md.log", "json-events") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_2" >> {
    "[aminoacids] test with json" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/aminoacids/md.log", "json") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_3" >> {
    "[argon] test with json-events" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/argon/md.log", "json-events") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_4" >> {
    "[argon] test with json" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/argon/md.log", "json") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_5" >> {
    "[water] test with json-events" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/water/md.log", "json-events") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_6" >> {
    "[water] test with json" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/water/md.log", "json") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_7" >> {
    "[tip4p] test with json-events" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/tip4p/md.log", "json-events") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_8" >> {
    "[tip4p] test with json" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/tip4p/md.log", "json") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_9" >> {
    "[Fe] test with json-events" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/fe_test/md.log", "json-events") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_10" >> {
    "[Fe] test with json" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/fe_test/md.log", "json") must_== ParseResult.ParseSuccess
    }
  }
}
