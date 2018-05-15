package eu.nomad_lab.parsers

import eu.nomad_lab
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import eu.{ nomad_lab => lab }
import scala.collection.breakOut

object GromacsParser extends SimpleExternalParserGenerator(
  name = "GromacsParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("GromacsParser")) ::
      ("parserId" -> jn.JString("GromacsParser" + lab.GromacsVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.GromacsVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """.*?\bLog\s*file\b.*\n(?:.*\n)*\s*:-\)\s*GROMACS\s*-\s*gmx\b.*""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/gromacs/parser/parser-gromacs/GromacsParser.py",
    "${mainFilePath}"),
  resList = Seq(
    "parser-gromacs/GromacsParser.py",
    "parser-gromacs/GromacsDictionary.py",
    "parser-gromacs/GromacsCommon.py",
    "parser-gromacs/GromacsCmdLineArgs.py",
    "parser-gromacs/setup_paths.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/gromacs.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-gromacs" -> "parsers/gromacs/parser/parser-gromacs",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info",
    "python" -> "python-common/common/python/nomadcore"
  ) ++ DefaultPythonInterpreter.commonDirMapping(),
  metaInfoEnv = Some(lab.meta.KnownMetaInfoEnvs.gromacs)
)

