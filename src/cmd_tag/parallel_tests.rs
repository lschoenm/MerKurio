use super::*;
use clap::Parser;

fn run_tag(
    input: &std::path::Path,
    output: &std::path::Path,
    threads: usize,
    extra: &[&str],
) -> Result<()> {
    let mut argv = vec![
        "merkurio".to_string(),
        "tag".into(),
        "-i".into(),
        input.display().to_string(),
        "-s".into(),
        "CTC".into(),
        "GAG".into(),
        "--threads".into(),
        threads.to_string(),
        "--chunk-size".into(),
        "7".into(),
        "-l".into(),
        output.with_extension("log").display().to_string(),
        "-j".into(),
        output.with_extension("json").display().to_string(),
    ];
    if !extra.contains(&"-S") {
        argv.extend(["-o".into(), output.display().to_string()]);
    }
    argv.extend(extra.iter().map(|s| s.to_string()));
    let crate::Commands::Tag(args) = crate::Cli::try_parse_from(argv)?.cmd else {
        unreachable!()
    };
    tag_records(args)
}

fn decoded(path: &std::path::Path) -> Result<Vec<u8>> {
    let (mut reader, header): (Box<dyn RecordReader>, bam::Header) =
        if path.extension().unwrap() == "bam" {
            let reader = bam::BamReader::from_path(path, 0)?;
            let header = reader.header().clone();
            (Box::new(reader), header)
        } else {
            let reader = bam::SamReader::from_path(path)?;
            let header = reader.header().clone();
            (Box::new(reader), header)
        };
    let mut bytes = Vec::new();
    {
        let mut writer = bam::sam::SamWriterBuilder::new().from_stream(&mut bytes, header)?;
        let mut record = bam::Record::new();
        while reader.read_into(&mut record)? {
            writer.write(&record)?;
        }
        writer.finish()?;
    }
    Ok(bytes)
}

#[test]
fn parallel_tag_matches_serial_across_formats_and_options() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let sam = dir.path().join("input.sam");
    let mut data = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:1\tLN:100000\n".to_string();
    for i in 0..150 {
        let seq = if i % 3 == 0 { "CTCCTCGAG" } else { "AAAAAAAAA" };
        data.push_str(&format!(
            "read{i}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\tFFFFFFFFF\tkm:Z:CTC,OLD\n"
        ));
    }
    fs::write(&sam, data)?;
    let bam = dir.path().join("input.bam");
    {
        let mut reader = bam::SamReader::from_path(&sam)?;
        let mut writer = bam::bam_writer::BamWriterBuilder::new()
            .additional_threads(0)
            .from_path(&bam, reader.header().clone())?;
        let mut record = bam::Record::new();
        while reader.read_into(&mut record)? {
            writer.write(&record)?;
        }
        writer.finish()?;
    }
    let options: &[&[&str]] = &[
        &[],
        &["-m"],
        &["-v"],
        &["--hash"],
        &["-a"],
        &["-q", "2"],
        &["-I"],
        &["-r"],
        &["-S"],
    ];
    for input in [&sam, &bam] {
        for extension in ["sam", "bam"] {
            for extra in options {
                let serial = dir.path().join(format!("serial.{extension}"));
                run_tag(input, &serial, 1, extra)?;
                for threads in [2, 3, 4] {
                    let parallel = dir.path().join(format!("parallel.{extension}"));
                    run_tag(input, &parallel, threads, extra)?;
                    if !extra.contains(&"-S") {
                        assert_eq!(decoded(&serial)?, decoded(&parallel)?);
                    }
                    let clean_log = |path: PathBuf| -> Result<Vec<String>> {
                        Ok(fs::read_to_string(path)?
                            .lines()
                            .skip(4)
                            .map(String::from)
                            .collect())
                    };
                    assert_eq!(
                        clean_log(serial.with_extension("log"))?,
                        clean_log(parallel.with_extension("log"))?
                    );
                    let clean_json = |path: PathBuf| -> Result<serde_json::Value> {
                        let mut value: serde_json::Value =
                            serde_json::from_str(&fs::read_to_string(path)?)?;
                        value.as_object_mut().unwrap().remove("meta_information");
                        Ok(value)
                    };
                    assert_eq!(
                        clean_json(serial.with_extension("json"))?,
                        clean_json(parallel.with_extension("json"))?
                    );
                }
            }
        }
    }
    Ok(())
}

#[test]
fn parallel_tag_empty_input_and_empty_filtered_batches() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("input.sam");
    for data in [
        "@HD\tVN:1.6\n".to_string(),
        format!(
            "@HD\tVN:1.6\n{}",
            "r\t4\t*\t0\t0\t*\t*\t0\t0\tAAA\tFFF\n".repeat(40)
        ),
    ] {
        fs::write(&input, data)?;
        for threads in [1, 2, 4] {
            let output = dir.path().join("out.sam");
            run_tag(&input, &output, threads, &["-m"])?;
            assert!(
                fs::read_to_string(output)?
                    .lines()
                    .all(|line| line.starts_with('@'))
            );
        }
    }
    Ok(())
}

#[test]
fn parallel_tag_propagates_invalid_tag_error() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("input.sam");
    fs::write(
        &input,
        format!(
            "@HD\tVN:1.6\n{}",
            "r\t4\t*\t0\t0\t*\t*\t0\t0\tCTC\tFFF\tkm:i:3\n".repeat(300)
        ),
    )?;
    for threads in [1, 2, 4] {
        let error = run_tag(&input, &dir.path().join("out.sam"), threads, &[]).unwrap_err();
        assert!(format!("{error:#}").contains("Invalid tag value format"));
    }
    Ok(())
}

#[test]
fn parallel_tag_without_logs_matches_serial() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("input.sam");
    fs::write(
        &input,
        format!(
            "@HD\tVN:1.6\n{}",
            "r\t4\t*\t0\t0\t*\t*\t0\t0\tCTCCTC\tFFFFFF\tkm:Z:OLD,CTC\n".repeat(80)
        ),
    )?;
    let mut expected = None;
    for threads in [1, 2, 3, 4, 0] {
        let output = dir.path().join("out.sam");
        let crate::Commands::Tag(args) = crate::Cli::try_parse_from([
            "merkurio",
            "tag",
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "-s",
            "CTC",
            "--threads",
            &threads.to_string(),
            "--chunk-size",
            "3",
        ])?
        .cmd
        else {
            unreachable!()
        };
        tag_records(args)?;
        let actual = fs::read_to_string(output)?;
        assert!(
            actual
                .lines()
                .filter(|line| !line.starts_with('@'))
                .all(|line| line.ends_with("km:Z:CTC,OLD"))
        );
        if let Some(expected) = &expected {
            assert_eq!(&actual, expected);
        } else {
            expected = Some(actual);
        }
    }
    Ok(())
}
