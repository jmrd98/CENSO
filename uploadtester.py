import os
import requests
import argparse
from urllib.parse import quote
from time import sleep

def upload_and_submit_files(
    files,
    upload_url="https://arrows.emsl.pnnl.gov//api/upload/",
    submit_url="https://arrows.emsl.pnnl.gov//api/submit_output",
):
    """
    Uploads multiple local files and submits them to a remote endpoint.

    Args:
        files (list): List of file paths to upload.
        upload_url (str): URL to upload the files.
        submit_url (str): URL to submit the uploaded files.

    Returns:
        str: Response message from the submission.
    """

    uploaded_files = []

    # Upload each file individually
    for file_path in files:
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue

        with open(file_path, "rb") as file:
            response = requests.post(upload_url, files={"file": file})
            sleep(1)  # Sleep for 1 second to avoid rate limiting

        if response.status_code == 200:
            print(f"Uploaded: {file_path}")
            uploaded_files.append(os.path.basename(file_path))
        else:
            print(f"Failed to upload: {file_path}, Status Code: {response.status_code}")

    if not uploaded_files:
        return "No files were uploaded."

    # URL-encode the filenames
    encoded_files = [quote(f) for f in uploaded_files]
    datafiles = " ".join(encoded_files)
    print(f"Submitting files: {datafiles}")
    # params = {"datafiles": datafiles}
    submit_response = requests.get(f"{submit_url}/{datafiles}")
    sleep(1)  # Sleep for 1 second to avoid rate limiting
    # submit_response = requests.get(submit_url, params=params)
    print(f"Submit URL: {submit_response.url}")
    print(f"Submit Response Status Code: {submit_response.status_code}")
    print(f"Submit Response Text: {submit_response.text}")

    if submit_response.status_code == 200:
        return submit_response.text
    else:
        return f"Failed to submit files, Status Code: {submit_response.status_code}"


if __name__ == "__main__":
    DEFAULT_UPLOAD_URL = "https://arrows.emsl.pnnl.gov/api/upload/"
    DEFAULT_SUBMIT_URL = "https://arrows.emsl.pnnl.gov/api/submit_output"

    parser = argparse.ArgumentParser(
        description="Upload and submit files to a remote endpoint."
    )
    parser.add_argument("--files", nargs="+", help="List of file paths to upload")
    parser.add_argument(
        "--upload_url",
        default=DEFAULT_UPLOAD_URL,
        help="URL to upload the files (default: %(default)s)",
    )
    parser.add_argument(
        "--submit_url",
        default=DEFAULT_SUBMIT_URL,
        help="URL to submit the uploaded files (default: %(default)s)",
    )

    args = parser.parse_args()

    response_message = upload_and_submit_files(
        args.files, args.upload_url, args.submit_url
    )
    print(response_message)
