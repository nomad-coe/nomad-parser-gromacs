/*
 * Copyright 2016-2018 Berk Onat, Fawzi Mohamed
 * 
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */

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
    "[tip4p] test with json-events" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/tip4p/md.log", "json-events") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_6" >> {
    "[tip4p] test with json" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/tip4p/md.log", "json") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_7" >> {
    "[Fe] test with json-events" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/fe_test/md.log", "json-events") must_== ParseResult.ParseSuccess
    }
  }
  "GromacsParserTest_8" >> {
    "[Fe] test with json" >> {
      ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/fe_test/md.log", "json") must_== ParseResult.ParseSuccess
    }
  }
  //"GromacsParserTest_9" >> {
  //  "[water] test with json-events" >> {
  //    ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/water/md.log", "json-events") must_== ParseResult.ParseSuccess
  //  }
  //}
  //"GromacsParserTest_10" >> {
  //  "[water] test with json" >> {
  //    ParserRun.parse(GromacsParser, "parsers/gromacs/test/examples/water/md.log", "json") must_== ParseResult.ParseSuccess
  //  }
  //}
}
